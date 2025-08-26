# test input

source("R/read_excel_allsheets.R")
library(readxl)
library("testthat")

test_that("read_excell_allsheets works if the file is the right format",{
  expect_error(read_excel_allsheets("Mouse_Test_Data.xlsx", tibble = FALSE), NA)
})

test_that("read_excell_allsheets detect group tab does not have colnames",{
  expect_error(object = read_excel_allsheets("inputs/Mouse_Test_Data_nocolnames.xlsx", tibble = FALSE), expected = "The last sheet of the excel file should contain a table with 2 columns and the same number of rows of the other sheets in the file")
})

test_that("read_excell_allsheets check that all sheet names are in the group table",{
  expect_error(object = read_excel_allsheets("inputs/Mouse_Test_Data different_names.xlsx", tibble = FALSE), expected = "all the experiment sheet names have to be reported in the group table")
})


test_that("read_excell_allsheets check that the string all is not contained in the sheet names",{
  expect_error(object = read_excel_allsheets("inputs/Mouse_Test_Data_all.xlsx", tibble = FALSE), expected = "the string all cannot be used in the sheet names")
})

test_that("read_excell_allsheets check that all the element of the second column of group table are numbers",{
  expect_error(object = read_excel_allsheets("inputs/Mouse_Test_Data_not_only_numbers.xlsx", tibble = FALSE), expected = "the second column of the group tab has to contain numbers")
})

test_that("read_excell_allsheets check that all the element of the second column of group table are consecutive numbers",{
  expect_error(object = read_excel_allsheets(filename = "inputs/Mouse_Test_Data consecutive_numbers.xlsx", tibble = FALSE), expected = "the numbers in the second column of the group tab should be consecutive")
})

test_that("read_excell_allsheets check the sheet are not empty",{
  expect_error(object = read_excel_allsheets(filename = "inputs/Mouse_Test_Data_empty.xlsx", tibble = FALSE), expected = "empty sheet")
})

test_that("read_excell_allsheets check the sheet are not empty",{
  expect_error(object = read_excel_allsheets(filename = "inputs/Mouse_Test_Data_empty2.xlsx", tibble = FALSE), expected = "empty sheet")
})
