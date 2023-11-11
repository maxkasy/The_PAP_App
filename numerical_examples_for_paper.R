library(kableExtra)
source("optimal_implementable_test.R")
source("optimal_simple_test.R")


binary_example = function(input){
    v_binary = list()
    v_binary$optimal_output = pap_binary_data(etaJ = input$etaJ,
                                              minp = input$minp_binary, size = input$size_binary,
                                              a = input$alpha, b = input$beta) 
    v_binary$optimal_output$test = v_binary$optimal_output$test |>
        filter(t>0) |>
        arrange(t)
    v_binary$simple_output = pap_binary_data_simple(etaJ = input$etaJ,
                                                    minp = input$minp_binary, size = input$size_binary,
                                                    a = input$alpha, b = input$beta)
    v_binary
}


normal_example = function(input){
    v_normal = list()
    v_normal$optimal_output = pap_normal_data(etaJ = input$etaJ, 
                                              mu0 = input$mu0, Sigma0 = input$Sigma0,
                                              mu = input$mu, Sigma = input$Sigma,
                                              size = input$size_normal)
    # v_normal$optimal_output$test = v_normal$optimal_output$test |>
    #     filter(t>0) |>
    #     arrange(t)
    v_normal$simple_output = pap_normal_data_simple(etaJ = input$etaJ, 
                                                    mu0 = input$mu0, Sigma0 = input$Sigma0,
                                                    mu = input$mu, Sigma = input$Sigma,
                                                    size = input$size_normal)
    v_normal
}     


binary_output = function(input_binary, v_binary, title){
    cat("\n\\paragraph{", title, "}\n")
    
    # cat('\n\nThe null hypothesis is that the mean is below ', input_binary$minp_binary, '. ',
    #     'The required size of the test is ', input_binary$size_binary, '. ', sep = "")
    
    cat('The probability of observing each of the components is (', paste0(input_binary$etaJ, collapse =", "), ').\\\\\n\n', sep = "")
    
    
    
    cat('\\begin{tabular}{c|c}\n',
        '\\textsc{Optimal test} &\\textsc{Optimal simple test}\\\\ \n',
        '\\toprule\n \\addlinespace\n',
        kable(v_binary$optimal_output$test, digits = 2,
              row.names = F, format = "latex", booktabs = TRUE),
        '&',
        kable(v_binary$simple_output$test, digits =2,
              row.names = F, format = "latex", booktabs = TRUE),
        '\\\\ \n  \\addlinespace\n',
        'Expected power: ',
        round(v_binary$optimal_output$expected_power, 3),
        '& Expected power: ',
        round(v_binary$simple_output$expected_power, 3),
        '\\end{tabular}', sep =""
    )
    
}

normal_output = function(input_normal, v_normal, title){
    figure_file = paste0("Figures/Normal-", title, ".png")
    ggsave(paste0("../../Manuscript/", figure_file),
           plot_2d_normal(v_normal$optimal_output$test),
           width = 4, height = 3)
    
    cat("\n\\paragraph{", title, "}\n")
    
    # cat('\n\nThe null hypothesis is that X has a mean vector of (', paste0(input_normal$mu0, collapse =", "), '). ',
    #     'The required size of the test is ', input_normal$size_normal, '. ', sep = "")
    
    cat('The probability of observing each of the components is (', paste0(input_normal$etaJ, collapse =", "), ').\n\n', sep = "")
    cat('\n\nThe interim prior is that X has a mean vector of (', paste0(input_normal$mu, collapse =", "), 
        '), and a variance of $\\bigl(\\begin{smallmatrix}', paste0(input_normal$Sigma[1,], collapse = " & "), '\\\\',
         paste0(input_normal$Sigma[2,], collapse = " & "), "\\end{smallmatrix}\\bigr)$.\\\\\n\n", sep = "")
    
    cat('\\begin{tabular}{c|c}\n',
    '\\textsc{Optimal test} &\\textsc{Optimal simple test}\\\\ \n',
    '\\toprule\n',
    '\\raisebox{-.5\\totalheight}{\\includegraphics[width=.5\\textwidth]{',
        figure_file,
    '}} &',
        kable(v_normal$simple_output$test, digits =2,
              row.names = F, format = "latex", booktabs = TRUE),
    '\\\\ \n',
    'Expected power: ',
        round(v_normal$optimal_output$expected_power, 3),
    '& Expected power: ',
        round(v_normal$simple_output$expected_power, 3),
    '\\end{tabular}', sep =""
    )

}




## Normal examples

file.remove("../../Manuscript/Sections/normal_examples.tex")
sink("../../Manuscript/Sections/normal_examples.tex")

input_normal = list(
    etaJ = c(.9,.5),
    mu0 = c(0,0),
    Sigma0 = matrix(c(1,0,0,1), nrow = 2),
    mu = c(1,1),
    Sigma = matrix(c(2,1,1,2), nrow = 2),
    size_normal = .05
) 
normal_output(input_normal, 
              normal_example(input_normal),
              "Example 1")

input_normal = list(
    etaJ = c(.9,.1),
    mu0 = c(0,0),
    Sigma0 = matrix(c(1,0,0,1), nrow = 2),
    mu = c(1,1),
    Sigma = matrix(c(2,1,1,2), nrow = 2),
    size_normal = .05
) 
normal_output(input_normal, 
              normal_example(input_normal),
              "Example 2")

input_normal = list(
    etaJ = c(.9,.9),
    mu0 = c(0,0),
    Sigma0 = matrix(c(1,0,0,1), nrow = 2),
    mu = c(1,1),
    Sigma = matrix(c(3,2,2,3), nrow = 2),
    size_normal = .05
) 
normal_output(input_normal, 
              normal_example(input_normal),
              "Example 3")


input_normal = list(
    etaJ = c(.9,.9),
    mu0 = c(0,0),
    Sigma0 = matrix(c(1,0,0,1), nrow = 2),
    mu = c(2,.5),
    Sigma = matrix(c(2,1,1,2), nrow = 2),
    size_normal = .05
) 
normal_output(input_normal, 
              normal_example(input_normal),
              "Example 4")

sink()




## Binary examples

file.remove("../../Manuscript/Sections/binary_examples.tex")
sink("../../Manuscript/Sections/binary_examples.tex")

input_binary = list(
    etaJ = c(.9,.5),
    minp_binary = .1,
    size_binary = .05,
    alpha = 1,
    beta = 1
)
binary_output(input_binary, 
              binary_example(input_binary),
              "Example 5")

input_binary = list(
    etaJ = c(.9,.5, .1),
    minp_binary = .1,
    size_binary = .05,
    alpha = 1,
    beta = 1
)
binary_output(input_binary, 
              binary_example(input_binary),
              "Example 6")


input_binary = list(
    etaJ = c(.9,.8,.7,.6),
    minp_binary = .1,
    size_binary = .05,
    alpha = 1,
    beta = 1
)
binary_output(input_binary, 
              binary_example(input_binary),
              "Example 7")

sink()


