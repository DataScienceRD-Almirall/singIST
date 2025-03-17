## code to prepare `example_superpathway_fit_model` object goes here
data <- example_superpathway_input
example_superpathway_fit_model <- fitOptimal(example_superpathway_input)

usethis::use_data(example_superpathway_fit_model, overwrite = TRUE)
