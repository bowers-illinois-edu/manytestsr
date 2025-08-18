# Functions for combining splitting and testing

#' Test, Split, Repeat
#'
#' Split and test.
#' @param idat Data at the unit level.
#' @param bdat Data at the block level.
#' @param blockid A character name of the column in idat and bdat indicating the block.
#' @param splitfn A function to split the data into two pieces --- using bdat
#' @param pfn A function to produce pvalues --- using idat.

#' @param alphafn A function to adjust alpha at each step. Takes one or more
#' p-values plus a stratum or batch indicator. Currently alpha_investing,
#' alpha_saffron, alpha_addis are accepted. All of them wrap the corresponding
#' functions from the `onlineFDR` package.

#' @param local_adj_p_fn Function. A function that adjusts p-values at a node (e.g. \code{local_simes}).

#' @param simthresh Below which number of total observations should the p-value
#' functions use permutations rather than asymptotic approximations

#' @param sims Number of permutations for permutation-based testing
#' @param maxtest Maximum splits or tests to do. Should probably not be smaller than the number of experimental blocks.

#' @param copydts TRUE or FALSE. TRUE if using find_blocks standalone. FALSE if
#' copied objects are being sent to find_blocks from other functions.

#' @param splitby A string indicating which column in bdat contains a variable
#' to guide splitting (for example, a column with block sizes or block harmonic
#' mean weights or a column with a covariate (or a function of covariates) or a
#' column with a factor with levels separated by "." that indicates a
#' pre-specified series of splits (see splitSpecifiedFactor))

#' @param stop_splitby_constant TRUE if the splitting should stop when splitby
#' is constant within a given branch of the tree. FALSE if splitting should
#' continue even when splitby is constant. Default is TRUE. Different
#' combinations of splitby, splitfn, and stop_splitby_constant make more or
#' less sense as described below.

#' @param blocksize A string with the name of the column in bdat contains
#' information about the size of the block (or other determinant of the power
#' of tests within that block, such as harmonic mean weight of the block or
#' variance of the outcome within the block.)

#' @param thealpha Is the error rate for a given test (for cases where alphafn is NULL, or the starting alpha for alphafn not null)

#' @param thew0 Is the starting "wealth" of the alpha investing procedure (this
#' is only relevant when alphafn is not null).

#' @param fmla A formula with outcome~treatment assignment  | block where
#' treatment assignment and block must be factors.

#' @param return_what Character. Return a data.table of blocks "blocks", a
#' data.table of nodes "nodes", or default both c("blocks","nodes").

#' @param final_global_adj Character. One of \code{"none"}, \code{"fdr"}, \code{"fwer"}.

#' @param parallel Should the pfn use multicore processing for permutation
#' based testing. Default is no. But could be "snow" or "multicore" following
#' `approximate` in the coin package.

#' @param ncores The number of cores used for parallel processing
#' @param trace Logical, FALSE (default) to not print split number. TRUE prints the split number.
#' @return A data.table containing information about the sequence of splitting and testing

#' @details Some notes about the splitting functions and how they relate to
#' splitting criteria (splitby) and stopping criteria (stop_splitby_constant).

#'
#'  * [splitCluster()] splits the blocks into groups that are as similar as
#' possible to each other on splitby using the kmeans clustering algorithm
#' (using a combination of [kmeans()] or [KMeans_rcpp()]). This will not work
#' with factor variables. When the splitting criteria is constant, it will
#' return random splits into roughly two equal sized groups of blocks if
#' stop_splitby_constant=FALSE. If stop_splitby_constant=TRUE then
#' [find_blocks()] will stop and return groups of blocks as detected or not.
#'  * [splitSpecifiedFactor()] will split the blocks into two groups following
#'  prespecified pattern encoded into the labels for the levels of the factor.
#'  For example, if we imagine three nested levels of splitting (like states,
#'  districts, neighborhoods), the factor would have labels like
#'  `category1_level1.category2_level1.category3_level1` and where splits will
#'  occur from left to right depending on whether there is existing variation
#'  at that level. When the factor is constant and stop_splitby_constant=TRUE
#'  splitting stops. For this reason we recommend that the right-most label of
#'  this factor be the individual blocks themselves---to ensure that testing
#'  descends to the block level if it can. When stop_splitby_constant=FALSE,
#'  then it uses random splits.
#'  * [splitSpecifiedFactorMulti()] will split the blocks into two or more
#'  groups following prespecified pattern encoded into the labels for the
#'  levels of the factor. For example, if we imagine three nested levels of
#'  splitting (like states, districts, neighborhoods), the factor would have
#'  labels like `category1_level1.category2_level1.category3_level1` and where
#'  splits will occur from left to right depending on whether there is existing
#'  variation at that level. For this reason we recommend that the right-most
#'  label of this factor be the individual blocks themselves---to ensure that
#'  testing descends to the block level if it can. When the factor is constant
#'  and stop_splitby_constant=TRUE splitting stops. When
#'  stop_splitby_constant=FALSE, then it uses random splits.
#'  * [splitEqualApprox()] splits the sets of blocks into two groups where the
#'  sum of the splitby vector is approximately the same in each split. For
#'  example, if splitby is number of units in a block, then this splitting
#'  function makes two groups of blocks, each group having the same total
#'  number of units. This splitting function will work with discrete or factors
#'  but will do: `rank_splitby <- rank(splitby)` and then divide the blocks
#'  into groups based on taking every other rank. So, for factors variables
#'  with few categories that are ordered, this will allocate every other
#'  category to one or another group.
#'  * [splitLOO()] chooses the blocks largest on the splitby vector one at a
#' time so that we have two tests, one focusing on the highest ranked block and
#' one on all of the rest of the blocks (for example, the block with the most
#' units in it versus the rest of the blocks). When the splitby vector has
#' ties, it chooses one block at random among those tied for the first or
#' largest rank. When the split vector has few values, for example, only two
#' values, it will still split assuming that the vector is numeric (so, 1 is
#' ranked higher than 0) and then randomly among ties. If
#' stop_splitby_constant=TRUE, then the algorithm will stop after exhausting
#' the blocks in the higher ranked category (thinking about the binary splitby
#' case). For this reason we advise against using splitLOO with a factor
#' splitby vector with few categories. [splitLOO()] is best used with a splitby
#' vector like block-size --- which could be constant and thus just create a
#' random choice of a single block or could vary and thus focus the testing on
#' the largest/highest ranked blocks.
#'
#' @importFrom stringi stri_count_fixed stri_split_fixed stri_split stri_sub stri_replace_all stri_extract_last
#' @importFrom digest digest getVDigest
#' @export
find_blocks <-
  function(idat,
           bdat,
           blockid = "block",
           splitfn,
           pfn,
           alphafn = NULL,
           local_adj_p_fn = NULL,
           simthresh = 20,
           sims = 1000,
           maxtest = 2000,
           thealpha = 0.05,
           thew0 = .05 - .001,
           fmla = YContNorm ~ trtF | blockF,
           parallel = "multicore",
           ncores = 4,
           copydts = FALSE,
           splitby = "hwt",
           stop_splitby_constant = TRUE,
           blocksize = "hwt",
           return_what = c("blocks", "nodes"),
           trace = FALSE) {
    # Some checks: We can split on a constant variable if we are collecting the
    # blocks into groups of equal size

    splitfn_text <- deparse(splitfn)
    split_fn_equal_approx <- length(grep("%%2", splitfn_text)) > 0
    split_fn_cluster <- length(grep("Ckmeans", splitfn_text) > 0)
    if (stop_splitby_constant & !split_fn_equal_approx) {
      stopifnot(
        "The splitby variable must have at least two values if stop_splitby_constant is TRUE" = uniqueN(bdat[[splitby]]) >= 2
      )
    }

    # Setup
    ## Initialize numeric node tracking
    node_tracker <- create_node_tracker()
    
    ## Should we always copy the data tables? And rename them?
    if (copydts) {
      bdat <- data.table::copy(bdat)
      idat <- data.table::copy(idat)
    } else {
      ## This is just some defensive programming to delete columns created in a previous run
      cols_to_del_b <-
        grep(
          "^p[0-9]|^g[0-9]|^alpha[0-9]|pfinal|testable|nodesize|blocksize|blocksbygroup|^nodenum|group_id",
          names(bdat),
          value = TRUE
        )
      cols_to_del_i <-
        grep(
          "^p[0-9]|^g[0-9]|^alpha[0-9]|pfinal|testable|nodesize|blocksize|blocksbygroup|^nodenum|group_id",
          names(idat),
          value = TRUE
        )
      if (length(cols_to_del_b) != 0) {
        bdat[, (cols_to_del_b) := NULL]
      }
      if (length(cols_to_del_i) != 0) {
        idat[, (cols_to_del_i) := NULL]
      }
    }
    # Some of the testing below needs the block id to be a factor
    # Step 1: Root of tree. One test. We are numbering simply following literature on k-ary trees where the root node is node 1
    i <- 1L
    bdat[, p1 := pfn(
      fmla = fmla,
      dat = idat,
      simthresh = simthresh,
      sims = sims,
      parallel = parallel,
      ncpu = ncores
    )]
    ## if(any(bdat$p1 < .05)){ browser() }   ## to look at false rejection processes
    bdat[, pfinalb := p1] # just to initialize pfinalb
    # Record the names of the group and pvalues.
    gnm <- "g1"
    pnm <- "p1"
    bdat[, group_id := 1L] # initialize with numeric group ID
    bdat[, biggrp := factor(rep("1", nrow(bdat)))] # for backward compatibility
    bdat[, g1 := biggrp]
    if (is.null(alphafn)) {
      bdat[, alpha1 := thealpha]
    } else {
      # Current thought: the FDR controlling procedure should not be more conservative
      # than the FWER controlling procedure.
      # So, use which ever alpha is bigger at the overall root.
      bdat[, alpha1 := max(
        alphafn(
          pval = unique(p1),
          batch = unique(g1),
          nodesize = sum(get(blocksize)),
          thealpha = thealpha,
          thew0 = thew0
        ),
        thealpha
      )]
    }
    # If p1 > alpha1 then stop testing and return the data.
    bdat[, testable := (pfinalb < alpha1)]
    # node_dat is a node level dataset
    node_dat <- data.table(
      parent = 0L,
      p = unique(bdat$p1),
      group_id = 1L,
      biggrp = factor("1"),
      a = unique(bdat$alpha1),
      batch = "p1",
      testable = unique(bdat$testable),
      nodenum = 1L,
      depth = 1L,
      nodesize = sum(bdat[, get(blocksize)])
    )
    bdat[, nodenum_current := 1L]
    bdat[, nodenum_prev := 0L]
    setkeyv(idat, blockid)
    setkeyv(bdat, "testable")
    bdat[, blocksbygroup := length(unique(get(blockid))), by = group_id]
    if (!all(bdat$testable)) {
      ## Return the objects
      if (all(return_what == "blocks")) {
        return(bdat)
      }
      if (all(return_what == "nodes")) {
        return(node_dat)
      }
      if (all(return_what %in% c("blocks", "nodes"))) {
        return(list(bdat = bdat, node_dat = node_dat))
      }
    }
    # Step 2: Iterate through the tree
    # Keep testing until no testing possible: criteria are (1) all blocks all size 1 and/or (2)
    # all p_i > alpha_i  OR (3) simulation
    # limits reached (for testing of the algorithm).
    # We might want at least 2 testable blocks.
    while (sum(bdat$testable, na.rm = TRUE) > 1 && i < maxtest) {
      # if(i==9){ browser() } ## for debugging
      i <- i + 1L
      if (trace) {
        message("Split number: ", i)
      }
      gnm <- paste0("g", i) # name of the grouping variable for the current split
      pnm <- paste0("p", i) # name of the p-value variable for the current split
      alphanm <- paste0("alpha", i)
      # Set Group to NA for blocks where we have to stop testing
      bdat[(testable), (gnm) := splitfn(bid = get(blockid), x = get(splitby)), by = group_id]
      
      # Update node tracking with numeric IDs
      bdat[(testable), nodenum_prev := nodenum_current]
      
      # Create unique node IDs for each unique split group value within each parent group
      testable_blocks <- bdat[(testable)]
      if (nrow(testable_blocks) > 0) {
        # Process each parent group separately to maintain proper tree structure
        parent_groups <- unique(testable_blocks$nodenum_current)
        
        for (parent_id in parent_groups) {
          parent_blocks <- testable_blocks[nodenum_current == parent_id]
          unique_split_values <- unique(parent_blocks[[gnm]])
          
          # Assign new node IDs for children of this parent
          for (k in seq_along(unique_split_values)) {
            split_val <- unique_split_values[k]
            new_node_id <- node_tracker$next_id
            node_tracker$next_id <- node_tracker$next_id + 1L
            
            # All blocks with this split value get the same node ID
            bdat[(testable) & nodenum_current == parent_id & get(gnm) == split_val, 
                 nodenum_current := new_node_id]
            
            # Add to tracker
            new_node <- data.table(
              node_id = new_node_id,
              parent_id = parent_id,
              depth = i
            )
            node_tracker$tracker <- rbindlist(list(node_tracker$tracker, new_node))
          }
        }
      }
      
      # Update group_id and create biggrp  
      bdat[(testable), group_id := nodenum_current]
      bdat[, biggrp := factor(group_id)]
      
      bdat[, nodesize := sum(get(blocksize)), by = group_id]
      # Now merge idat and bdat again since the test has to be at the idat level
      idat[bdat, c("testable", "group_id", "biggrp") := mget(c("i.testable", "i.group_id", "i.biggrp")), on = blockid]
      pb <- idat[(testable), list(
        p = pfn(
          fmla = fmla,
          dat = .SD,
          simthresh = simthresh,
          sims = sims,
          parallel = parallel,
          ncpu = ncores
        )
      ), by = group_id]
      pb[, depth := i]
      # Use numeric node IDs directly - no string splitting needed
      pb[, nodenum := group_id]
      pb[, parent := get_parent_from_tracker(node_tracker, group_id)]
      # call "blocksize" the sum of the block sizes within group
      pb[bdat, nodesize := i.nodesize, on = "group_id"]
      ## Adjust the p-values for a given node
      if (!is.null(local_adj_p_fn)) {
        pb[, p := local_adj_p_fn(p), by = parent]
      }
      # Add biggrp to pb for compatibility
      pb[, biggrp := factor(group_id)]
      
      node_dat <-
        rbind(node_dat[, .(parent, p, a, group_id, biggrp, batch, testable, nodenum, depth, nodesize)],
          pb[, batch := pnm],
          fill = TRUE
        )
      setnames(pb, "p", pnm)
      setkeyv(pb, "group_id")
      bdat[pb, (pnm) := get(paste0("i.", pnm)), on = "group_id"]
      # bdat[(testable), pfinalb := pmax(get(pnm), pfinalb)]
      bdat[(testable), pfinalb := get(pnm)]
      # Now decide which blocks (units) can be tested again.
      # If a split contains only one block. We cannot test further.
      bdat[, blocksbygroup := .N, by = group_id]
      if (is.null(alphafn)) {
        bdat[, (alphanm) := thealpha]
        # Recall that `:=` **updates** values. So, it doesn't overwrite (testable==FALSE) values
        bdat[, testable := fifelse((pfinalb <= get(alphanm)) & (blocksbygroup > 1), TRUE, FALSE)]
        node_dat[, (alphanm) := thealpha]
      } else {
        if (i == 2) {
          node_dat[, (alphanm) := alphafn(
            pval = p,
            batch = batch,
            nodesize = nodesize,
            thealpha = thealpha,
            thew0 = thew0
          )]
          node_dat[is.na(a), a := get(alphanm)]
        } else {
          ## Find the nodes that are ancestors and descedents of each other since
          ## the alpha adjusting depends on this order

          ## mid_roots <- node_dat[depth == (i - 1) & (testable), nodenum]

          find_paths <- function(j) {
            tmp <- node_dat[depth == j, ]
            # Build paths using numeric ancestry from tracker
            paths <- lapply(tmp$nodenum, function(node_id) {
              build_numeric_ancestry(node_tracker, node_id)
            })
            return(paths)
          }

          setkey(node_dat, nodenum)
          thepaths <- find_paths(j = i)
          ## Do alpha adjustment for each path through the binary tree
          for (j in seq_along(thepaths)) {
            node_dat[J(thepaths[[j]]), (alphanm) := alphafn(
              pval = p,
              batch = depth,
              nodesize = nodesize,
              thealpha = thealpha,
              thew0 = thew0
            )]
          }
        }

        node_dat[is.na(a), a := get(alphanm)]
        setkey(node_dat, group_id)
        bdat[node_dat, (alphanm) := get(paste0("i.", alphanm)), on = "group_id"]
        # In deciding which blocks can be included in more testing we use
        # current p rather than max of previous p. A block is testable if current
        # p <= the alpha level AND the number of blocks in the group containing
        # the block is more than 1 We stop testing when we have only a single
        # block within a group (or branch) because the block is the unit.

        bdat[, testable := fifelse((get(pnm) <= get(alphanm)) &
          (blocksbygroup > 1), TRUE, FALSE)]

        # Also stop testing for groups within which we cannot split any more for
        # certain splitters. Currently set by hand.
      }

      if (stop_splitby_constant || split_fn_cluster) {
        ## Here there is no sense in spliting by differences in covariate value
        ## (creating clusters using k-means) if covariate values do not differ

        bdat[, testable := fifelse(uniqueN(get(splitby)) == 1, FALSE, unique(testable)), by = group_id]
      }
      node_dat[is.na(a), a := get(alphanm)]
      setkeyv(bdat, "testable") # for binary search speed
      # message(paste(unique(signif(bdat$pfinalb,4)),collapse=' '),appendLF = TRUE)
      # message('Number of blocks left to test: ', sum(bdat$testable))
    }
    ## Return the objects
    if (all(return_what == "blocks")) {
      return(bdat)
    }
    if (all(return_what == "nodes")) {
      return(node_dat)
    }
    if (all(return_what %in% c("blocks", "nodes"))) {
      return(list(bdat = bdat, node_dat = node_dat))
    }
  }


# Numeric node tracking system for performance

#' Create node tracker object
#' @return List with tracker data.table and next_id counter
#' @keywords internal
create_node_tracker <- function() {
  list(
    tracker = data.table(
      node_id = 1L,
      parent_id = 0L,
      depth = 1L
    ),
    next_id = 2L
  )
}

#' Add nodes to tracker
#' @param tracker Node tracking object
#' @param parent_ids Vector of parent node IDs
#' @param depth Current depth
#' @param n_children Number of children to create
#' @return Updated tracker with new node IDs
#' @keywords internal
add_nodes_to_tracker <- function(tracker, parent_ids, depth, n_children) {
  new_ids <- seq.int(tracker$next_id, length.out = n_children)
  tracker$next_id <- tracker$next_id + n_children
  
  new_nodes <- data.table(
    node_id = new_ids,
    parent_id = rep(parent_ids, length.out = n_children),
    depth = depth
  )
  
  tracker$tracker <- rbindlist(list(tracker$tracker, new_nodes))
  list(tracker = tracker, new_ids = new_ids)
}

#' Get parent ID from tracker
#' @param tracker Node tracking object  
#' @param node_ids Vector of node IDs to find parents for
#' @return Vector of parent node IDs (0 for root)
#' @keywords internal
get_parent_from_tracker <- function(tracker, node_ids) {
  sapply(node_ids, function(nid) {
    parent <- tracker$tracker[node_id == nid, parent_id]
    if (length(parent) == 0) 0L else parent
  })
}

#' Build ancestry path for numeric node ID
#' @param tracker Node tracking object
#' @param node_id Node ID to build path for
#' @return Integer vector of ancestry path from root to node
#' @keywords internal
build_numeric_ancestry <- function(tracker, node_id) {
  path <- c()
  current <- as.integer(node_id)
  while (current != 0) {
    path <- c(current, path)
    current <- get_parent_from_tracker(tracker, current)
  }
  return(path)
}

#' Use hashing to make a node id (DEPRECATED)
#'
#' This function is deprecated in favor of numeric node tracking.
#' Kept for backward compatibility only.
#'
#' @param d A vector
#' @return a vector of hashes
#' @importFrom digest digest getVDigest
#' @export
nodeidfn <- function(d) {
  crc32hash <- getVDigest(algo = "crc32")
  crc32hash(d)
}
