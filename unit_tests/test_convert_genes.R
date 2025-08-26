# test convert_genes function

source("R/convert_genes.R")
source("R/read_excel_allsheets.R")

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(AnnotationDbi)
library("testthat")


inputs = read_excel_allsheets("Mouse_Test_Data.xlsx", tibble = FALSE)
GList = inputs[[1]]
test_that("conveert_genes works fine",{
  expect_error(convert_genes(GList, organism = "Mouse", annType = "SYMBOL"), NA)
})

test_that("conveert_genes wrong organism",{
  expect_error(convert_genes(GList, organism = "organism", annType = "SYMBOL"), "Organism value must be one between Human, Mouse or Rat!")
})

test_that("conveert_genes wrong annoType",{
  expect_error(convert_genes(GList, organism = "Rat", annType = "type"), "annType value must be one between SYMBOL, ENSEMBL or ENTREZID!")
})

test_that("conveert_genes. wrong organism, the input genes comes from Rat",{
  expect_error(convert_genes(GList, organism = "Human", annType = "SYMBOL"), "the genes do not maps the organism or you are specifying the wrong annotation type")
})

test_that("conveert_genes wrong annotation type, the input file contains SYMBOLS",{
  expect_error(convert_genes(GList, organism = "Rat", annType = "ENTREZID"), "the genes do not maps the organism or you are specifying the wrong annotation type")
})

test_that("conveert_genes wrong annotation type, the input file contains SYMBOLS",{
  expect_error(convert_genes(GList, organism = "Rat", annType = "ENSEMBL"), "the genes do not maps the organism or you are specifying the wrong annotation type")
})
