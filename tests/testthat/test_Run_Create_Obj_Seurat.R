test_that("Run_Create_Obj returns a valid Seurat object for Slide-seq", {
  
  gene.count <- matrix(rpois(5 * 3, lambda = 20), nrow = 5)
  rownames(gene.count) <- paste0("Gene", seq_len(5))
  colnames(gene.count) <- paste0("Spot", seq_len(3))

  matched.data <- data.frame(
    spatial_name    = colnames(gene.count),
    barcode_sequence = paste0("BC", seq_len(3)),
    X_coordinate     = runif(3, 0, 100),
    Y_coordinate     = runif(3, 0, 100),
    UMI_count        = sample(50:200, 3),
    stringsAsFactors = FALSE
  )

  seu <- Run_Create_Obj(
    gene.matrix = gene.count,
    matched.data = matched.data,
    obj.type    = "Seurat",
    tech        = "Slideseq"
  )
  
  # check object class
  expect_s4_class(seu, "Seurat")
  expect_true("Spatial" %in% Seurat::Assays(seu))

  # check metadata 
  md <- seu@meta.data
  expect_true(all(c("X", "Y", "UMI") %in% colnames(md)))
  expect_equal(rownames(md), matched.data$barcode_sequence)
  
  # check coordinate information
  imgs <- names(seu@images)
  coords <- if ("image" %in% imgs) seu@images$image@coordinates else seu@images$fov@coordinates
  expect_equal(rownames(coords), matched.data$barcode_sequence)
  expect_equal(as.numeric(coords[,"X"]), matched.data$X_coordinate)
  expect_equal(as.numeric(coords[,"Y"]), matched.data$Y_coordinate)

})
