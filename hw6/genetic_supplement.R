## Supplemental code for homework 6
# Author: Eric Kalosa-Kenyon

model1 = function(x){
    t = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
    h = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
    temp = 0
    for(i in 1:11) {
        temp = temp + h[i]/2 * (1 + sign(x - t[i]))
    }
    return(temp)
}

noise1 = function(x){
 return(rnorm(length(x))/3)
}

random_chromosome = function(num_genes, nucleobases=c(0,1), max_non0=30){
    # Generates a chromosome
    #   lenght of chromosome is (num_genes)
    #   each element of the chromosome is an element of (nucleobases)
    non0 = sample(1:num_genes, max_non0)
    x = rep(nucleobases[1], num_genes)
    x[non0] = sample(nucleobases, max_non0, replace=TRUE)
    x = paste(x, collapse='')
    return(x)
}

mate_chromosomes = function(parents, xr=0.1, mr=0.1){
    # Produce a single offspring from parent chromosomes (p1, p2) with a
    # crossover rate (xr) and a mutation rate (mr)

    p1 = parents[1]
    p2 = parents[2]
    stopifnot(nchar(p1) == nchar(p2))

    child = ""
    w = 1 # crossover parameter
    for(i in 1:nchar(p1)){

        # Give child parents' genes for each gene (i)
        #   from p1 if w == 1
        #   from p2 if w == 0
        if(w){
            chr = substr(p1, i, i) # get i'th character from parent 1
        }else{
            chr = substr(p2, i, i) # get i'th character from parent 2
        }

        # Mutate gene
        if(runif(1)<mr){
            chr = as.character((as.integer(chr) + 1) %% 2)
        }

        # Extend child with the i'th parental character
        child = paste(child, chr, sep="")

        # Crossover if required
        if(runif(1)<xr){
            w = (w + 1) %% 2
        }
    }
    return(child)
}

regression_loss = function(reg, raw){
    # Sum of square error betweeen regression and raw
    # len(reg) == len(raw)
    return(sum((reg-raw)^2))
}

regress_chromosome = function(chrom, raw){
    # return c(fhat(x))
    # chrom is a string of "0001110101.."
    # raw is a c(num), len(chrom) == len(raw) - 1

    # Initialize working objects
    segment_observations = c(raw[1])
    regression = c()

    for(i in 1:nchar(chrom)){ # for each gene in the chromosome
        g = as.integer(substr(chrom, i, i)) # extract the gene from the string

        if(g){ # if the gene indicates a breakpoint
            # calculate the regression for that segment

            regression = c(regression,
                           rep(mean(segment_observations),
                               length(segment_observations)
                               )
                           )
            segment_observations = c(raw[i+1]) # reset to next segment

        }else{ # if the gene does not indicate a breakpoint
            # extend the segment
            segment_observations = c(segment_observations, raw[i+1])
        }
    }

    # calculate final segment regression
    regression = c(regression,
                   rep(mean(segment_observations),
                       length(segment_observations)
                       )
                   )

    return(regression)
}

calc_AIC = function(n, k, loss){
    # AIC is k-ln(sum sq losses)
    #   n is positive integer
    #   k is num params
    #   loss is losses
    #   k and loss can be c() as long as they have the same length
    # Note that this AIC can be alternatively parameterized

    r = 2*k + log(loss)
    return(r)
}

count_params = function(vec_of_chromosomes){
    # given a vector of chromosomes (each is "00111010111...")
    #   count the number of segments this represents i.e. the number of
    #   parameters that are implicitly contained in this chromosome.

    r = sapply( # vector of integers summing number of breakpoints
            lapply(
                strsplit(vec_of_chromosomes, ""),
                as.integer),
            sum) + 1 # number of constants is breaks + 1
    return(r)
}

fitness = function(population, raw){
    # Calculate fitness for each individual in the population relative to the
    #   raw data.
    #   Complexity is either "AIC" or "MDL". Raw is the raw data to fit.
    # Return a vector of floats between 0 and 1.
    #   The AIC and MDL are normalized to [0,1]

    # Check that input is ok
    stopifnot(nchar(population$Chromosomes[1]) + 1 == length(raw))

    # Initialize parameters
    losses = c()
    K = dim(population)[1]  # (K) is number of individuals in population
    n = length(raw) # (n) is the number of data observations

    # Calculate losses for each individual in the generation
    for(k in 1:K){
        x = population$Chromosome[k] # extract a chromosome
        regr = regress_chromosome(x, raw) # calculate the regression
        loss = regression_loss(regr, raw) # and calculate the loss
        losses = c(losses, loss) # record that individual's loss
    }

    ## Penalize for complexity
    ns_params = count_params(population$Chromosome)
    AIC = calc_AIC(n, ns_params, losses)

    # calculate sum of log(\hat{n_j}) for MDL calculation
    # sum_n_hat = for each chromosome, length of 0s between 1s plus 1
    # TODO: refactor into a separate subroutine
    #   MDL = calc_MDL(n, ns_params, losses)
    sum_n_hat = sapply(
                    lapply( # for each chromosome, get \hat{n_j}'s
                        strsplit( # get strs of 0's (whose lengths = n_j-1)
                            population$Chromosome,
                            "1"
                        ),
                        nchar),
                    (function (x) sum(log(x+1))) # sum(log(\hat{n_j}))
                )
    MDL = ns_params*log(n) +
        1/2*sum_n_hat +
        n/2*log(losses/n)

    return(list(AIC, MDL))
}

plot_chromosome = function(chrom, raw_df, color="blue"){
    # return a ggplot object with piecewise constants plotted on the scatter

    fhat = regress_chromosome(chrom, raw_df$y)
    fhat_df = data.frame(x=raw_df$x, y=fhat)

    plt = ggplot(raw_df, aes(x=x, y=y)) +
        geom_point() +
        geom_line(data=fhat_df, color=color) +
        xlab("Index") + ylab("Value")

    return(plt)
}

plot_generation = function(gen, raw_df, wh_score="AIC"){
    # return a ggplot object with piecewise constants plotted on the scatter
    #   for each chromosome in the generation
    # Currently, colorized lines based on AIC/MDL are unimplemented

    wh_gen = gen$Generation[1]

    # prepare base plot
    plt = ggplot(raw_df, aes(x=x, y=y)) +
        geom_point() +
        xlab("Index") + ylab("Value") +
        ggtitle(paste("Generation", wh_gen, wh_score))

    # prepare regression data for plotting
    fhat_df = data.frame(matrix(ncol = 5, nrow = 0))
    colnames(fhat_df) = c("x","y","Group","AIC","MDL")
    for(i in 1:dim(gen)[1]){
        chrom = gen$Chromosome[i]
        fhat = regress_chromosome(chrom, raw_df$y)
        fhat_df_i = data.frame(x=raw_df$x, y=fhat)
        fhat_df_i$Group = i
        fhat_df_i$AIC = as.numeric(gen$AIC[i])
        fhat_df_i$MDL = as.numeric(gen$MDL[i])
        fhat_df = rbind(fhat_df, fhat_df_i)
    }

    # add geomlines for aic or mdl depending on which was passed into this fxn
    if(wh_score=="AIC"){
        plt = plt + geom_line(data=fhat_df, aes(color=AIC, group=Group))
    }else if(wh_score=="MDL"){
        plt = plt + geom_line(data=fhat_df, aes(color=MDL, group=Group))
    }

    # plot AIC or MDL in a color gradient
    plt = plt + scale_colour_gradientn(colours=rainbow(4))

    return(plt)
}

best_score_and_org = function(population, wh_score){
    # Determine best score in the population and which organism achieved it
    if(wh_score=="AIC"){
        best_score = min(as.numeric(population$AIC))
        which_organism_is_best = which(
            as.numeric(population$AIC) == best_score)[1]
    }else if(wh_score=="MDL"){
        best_score = min(as.numeric(population$MDL))
        which_organism_is_best = which(
            as.numeric(population$MDL) == best_score)[1]
    }

    return(list(best_score, which_organism_is_best))
}
