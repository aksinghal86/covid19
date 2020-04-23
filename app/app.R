#### Libraries ####
library(shiny)
library(shinythemes)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(scales)
library(gridExtra)
library(lubridate)
library(incidence)
library(epitrix)
library(distcrete)
library(EpiEstim)
library(earlyR)
library(minpack.lm)

#### Import data ####

# Public health interventions
interventions <- read.delim('data/interventions.csv', sep = ",", header = TRUE, stringsAsFactors = FALSE) %>% 
    mutate(stay_at_home = mdy(stay_at_home), 
           ed_facilities = mdy(ed_facilities), 
           non_essential = mdy(non_essential))

# COVID19 data
state_data <- read.delim('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv', sep=",", header = TRUE, stringsAsFactors = FALSE) %>%
#state_data <- read.delim('data/us-states.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(date = ymd(date), cum_cases = cases, cum_deaths = deaths) %>%
    group_by(state) %>% 
    mutate(cases = c(cum_cases[1], diff(cum_cases)), 
           deaths = c(cum_deaths[1], diff(cum_deaths)), 
           county = "All counties")

all_us <- state_data %>% 
    group_by(date) %>% 
    summarize(cases = sum(cases), deaths = sum(deaths), 
              cum_cases = sum(cum_cases), cum_deaths = sum(cum_deaths)) %>% 
    mutate(state = "All states", fips = 0, county = 'All counties')


county_data <- read.delim('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv', sep=",", header = TRUE, stringsAsFactors = FALSE) %>%
# county_data <- read.delim('data/us-counties.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(date = ymd(date), cum_cases = cases, cum_deaths = deaths) %>% 
    group_by(state, county) %>% 
    mutate(cases = c(cum_cases[1], diff(cum_cases)), 
           deaths = c(cum_deaths[1], diff(cum_deaths)))

all_data <- dplyr::bind_rows(all_us, state_data, county_data) %>% 
    left_join(interventions, by = "state") %>%
    filter(cases >= 0)

# AQI data
aqi_data <- read.delim('data/air_quality_data.csv', sep = ",", header = TRUE, stringsAsFactors = FALSE) %>%
    left_join(interventions, by = "state")

#### Parameters and functions for projections ####

# Set SI parameters for estimating R. Default from from Du et al. (2020)
si_distr_params <- list('mu' = 3.96, 'sigma' = 4.75)

# Create incidence objects that will be used for plotting 
create_incidence_objects <- function(dat, peak_date) {
    inc_cases <- dat %>%
        ungroup(state) %>%
        select(date, cases) %>%
        uncount(cases)

    # Create an incidence object
    i <- incidence(inc_cases$date)
    
    # If the user has provided a peak date, use that; otherwise
    # use the  peak date found by the model. The default peak date
    # is today + 1, so that the model can initialize and find the 
    # ideal peak date (a workaround, not an error!)
    if (peak_date == (today() + 1)) {
        i.peak <- incidence::find_peak(i)
    } else {
        i.peak <- peak_date
    }
    i.fit <- incidence::fit(subset(i, to = i.peak))
    
    return(list('i' = i, 
                'i.peak' = i.peak, 
                'i.fit' = i.fit))
}

# Simple logistic regression projections
tfwd = 60 # number of days in the future to forecast for

# Logistic regression model using nlsLm from minpack.lm package
simple_logistic_regression <- function(t, y) {
    fit <- nlsLM(y ~ p/(1 + exp(-a * (t-B))), start = list(p=1000, a=0.18, B=20), lower = c(min(y), 0, 0))
    return(fit)
}

# Function to make a prediction db from logistic regression fit
make_preds <- function(dat, type) {
    if (type == "cases") {
        dat <- dat %>% ungroup %>% select(date, cum_cases, cases)
        t <- 1:nrow(dat)
        y <- dat$cum_cases
    } else if (type == "deaths") {
        dat <- dat %>% ungroup %>% select(date, cum_deaths, deaths)
        t <- 1:nrow(dat)
        y <- dat$cum_deaths
    } else{
        return("Error! Incorrect type provided; type must be either 'cases' or 'deaths'")
    }
    
    # Fit regression model
    mod_fit <- simple_logistic_regression(t, y)
    
    # Make predctions and create a database
    Day <- 1:(nrow(dat) + tfwd)
    preds <- data.frame(date = dat$date[1] + days(Day-1)) %>%
        mutate(cum_pred = ceiling(predict(mod_fit, list(t = Day))),
               inc_pred = c(cum_pred[1], diff(cum_pred))) %>%
        left_join(dat, by = "date")

    return(list('mod_fit'=mod_fit, 'preds'=preds))
}

#### Plotting functions ####
# Set default theme 
my_theme <- theme_minimal() +
    theme(legend.position = 'none',
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 14),
          plot.subtitle = element_text(size = 12),
          panel.grid.major = element_line(color = "darkgrey")
    )

### Plots
# Create log cumulative plot
cum_log_plot <- function(dat, county, state) {
    ggplot(dat, aes(x = date)) +
        geom_line(aes(y = log(cum_cases)), color = "turquoise", size = 1) +
        geom_line(aes(y = log(cum_deaths)), color = "turquoise", size = 1, linetype = "dashed") +
        labs(x = "", y = "N (log scale)",
             title = paste0("Cumulative cases and deaths in ", county, ", ", state),
             subtitle = "(solid line = confirmed cases, dashed line = confirmed deaths)") +
        my_theme 
}

epi_plot <- function(dat, county, state) {
    ggplot(dat, aes(x = date)) +
        geom_line(aes(y = cases), color = "orange",  size = 1) +
        geom_line(aes(y = deaths*10), color = "orange", size = 1, linetype = "dashed") +
        # annotate(geom = "text", 
        #          x = date, y = max(dat$cases, na.omit = TRUE), 
        #          hjust = 0, vjust = 2,  
        #          color = 'red', size = 5, 
        #          label = 'Deaths MULITPLIED by 10\nfor visualization purposes\nNOT ACTUAL DEATHS!!!') +
        labs(x = "", y = "N", 
             title = paste0("Confirmed daily cases in ", county, ", ", state), 
             subtitle = "(solid line = confirmed cases, dashed line = confirmed deaths * 10 (SCALED!))") +
        my_theme + 
        scale_y_continuous(labels = comma_format(accuracy = 1))
}

# Create loglinear fit plot
loglinear_fit_plot <- function(inc_obj, county, state) {
    plot(inc_obj$i, color="turquoise", border="white") %>%
        add_incidence_fit(inc_obj$i.fit) +
        labs(title = paste0("Observed and modeled incidence in ", county, ", ", state)) +
        scale_x_date(date_labels = "%b %d") +
        my_theme + 
        scale_y_continuous(labels = comma_format())
}

# R0 distribution plot
r0_distr_plot <- function(inc_obj, county, state) {
    set.seed(1)
    
    # Set mu and sigma
    mu <- si_distr_params$mu
    sigma <- si_distr_params$sigma
    
    # Get the discrete gamma distribution
    param <- epitrix::gamma_mucv2shapescale(mu, sigma/mu)
    w <- distcrete("gamma", interval = 1, shape = param$shape, scale = param$scale, w = 0)
    
    # Sample growth R0 based on the SI distribution
    growth_R0 <- lm2R0_sample(inc_obj$i.fit$model, w)
    
    # plot
    hist(growth_R0, col = "turquoise", border = "white", main = "Distribution of growth R0", xlab = '') + 
        my_theme
}

# Instantaneous reproductive rate plot
re_plot <- function(dat, county, state) {
    
    # Set mu and sigma
    mean_si = si_distr_params$mu
    std_si = si_distr_params$sigma
    
    # Estimate R using the EpiEstim package
    res <- EpiEstim::estimate_R(
        dat %>% select(date, cases) %>% rename(dates = date, I = cases), 
        method = "parametric_si", config = make_config(list(mean_si = mean_si, std_si = std_si))
    )
    
    # Plot the instantaneous R
    plot(res, "R", xlab = '') +
        geom_hline(yintercept = 1, color = "red", size = 1) +
        geom_vline(xintercept = unique(dat$ed_facilities), lty = 'dashed', color = 'darkgreen', size = 1) +
        geom_vline(xintercept = unique(dat$non_essential), lty = 'dashed', color = 'darkblue', size = 1) +
        geom_vline(xintercept = unique(dat$stay_at_home), lty = 'dashed', color = 'darkorange', size = 1) +
        labs(x = "", y = "R", title = paste0("Estimated instantaneous R for ", county, ", ", state),
             subtitle = "(green dashed line = Educational services closed,\nblue dashed line = Non-essential services closed, \norange dashed line = Stay at home order)") +
        my_theme +
        theme(legend.position = 'none')
}

# Force of infection (lambda plot)
lambda_plot <- function(inc_obj, county, state) {
    
    # Set mu and sigma
    si_mean = si_distr_params$mu
    si_sd = si_distr_params$sigma
    simple_R <- earlyR::get_R(inc_obj$i, si_mean = si_mean, si_sd = si_sd, max_R = 8)
    plot(simple_R, 'lambdas', bty = 'n')
    abline(v = today(), col = 'darkorange', lty = 'dashed', lwd = 2)
}

# Projection plots for daily and cumulative cases using simple logistic regression
logistic_regression_fits <- function(dat, county, state) {
    cases_preds <- make_preds(dat, type = "cases")
    deaths_preds <- make_preds(dat, type = "deaths")

    # Cumulative plot for cases with projections
    cum_cases_plot <- ggplot(cases_preds$preds, aes(x=date)) +
        geom_line(aes(y=cum_pred), color = "turquoise", size = 2) +
        geom_point(aes(y=cum_cases), color = "darkcyan", shape = 1, size = 3) +
        geom_line(aes(y=cum_cases), color = "gray10", linetype = "dashed") +
        labs(x = "", y = 'Cumulative Incidence',
             title = paste('Cumulative cases fitted vs observed in', county, ",", state),
             subtitle = '(line = fitted incidence, circle = observed incidence)') +
        my_theme +
        scale_y_continuous(labels = comma_format())
    
    # Epi plot for cases with projections
    epi_cases_plot <- ggplot(cases_preds$preds, aes(x=date)) +
        geom_line(aes(y=inc_pred), color = "turquoise", size = 2) +
        geom_point(aes(y=cases), color = "darkcyan", shape = 1, size = 3) +
        geom_line(aes(y=cases), color = "gray10", linetype = "dashed") +
        labs(x = "", y = 'Daily incidence',
             title = paste('Daily cases fitted vs observed in', county, ",", state),
             subtitle = '(solid line = fitted incidence, circle = observed incidence)') +
        my_theme +
        scale_y_continuous(labels = comma_format())
    
    # Cumulative plot for deaths with projections
    cum_deaths_plot <- ggplot(deaths_preds$preds, aes(x=date)) +
        geom_line(aes(y=cum_pred), color = "orange", size = 2) +
        geom_point(aes(y=cum_deaths), color = "tan1", shape = 1, size = 3) +
        geom_line(aes(y=cum_deaths), color = "gray10", linetype = "dashed") +
        labs(x = "", y = 'Cumulative Incidence',
             title = paste('Cumulative deaths fitted vs observed in', county, ",", state),
             subtitle = '(solid line = fitted incidence, circle = observed incidence)') +
        my_theme +
        scale_y_continuous(labels = comma_format())
    
    # Epi plot for deaths with projections
    epi_deaths_plot <- ggplot(deaths_preds$preds, aes(x=date)) +
        geom_line(aes(y=inc_pred), color = "orange", size = 2) +
        geom_point(aes(y=deaths), color = "tan1", shape = 1, size = 3) +
        geom_line(aes(y=deaths), color = "gray10", linetype = "dashed") +
        labs(x = "", y = 'Daily incidence',
             title = paste('Daily deaths vs observed in', county, ",", state),
             subtitle = '(solid line = fitted incidence, circle = observed incidence)') +
        my_theme +
        scale_y_continuous(labels = comma_format())
    
    return(list('cases_fit' = cases_preds$mod_fit, 'deaths_fit' = deaths_preds$mod_fit,
                'cum_cases_plot' = cum_cases_plot, 'epi_cases_plot' = epi_cases_plot,
                'cum_deaths_plot' = cum_deaths_plot, 'epi_deaths_plot' = epi_deaths_plot))
}

aqi_plot <- function(dat, county, state, by = "day") {
    regions <- data.frame(ystart = c(0, 50), yend = c(50, 100), col = c("Good Air Quality", "Moderate Air Quality"))
    
    g <- ggplot() + 
        geom_rect(data = regions, aes(xmin = -Inf, xmax = Inf, ymin = ystart, ymax = yend, fill = col), alpha = 0.3) +
        scale_fill_manual(values = c("green1", "yellow1")) +
        labs(x = "", y = "AQI for PM2.5", fill = "") +
        theme(legend.position = 'bottom', legend.background = element_rect(fill = 'gray95')) 

    # Create plots by month, week or day based on user input
    if (by == "month") {
        subdat <- dat %>%
            mutate(month = factor(month.abb[month], levels = month.abb)) %>%
            group_by(year, month) %>%
            summarize(monthly_aqi = mean(daily_aqi), monthly_PM2.5 = mean(daily_PM2.5))
        
        g +
            new_scale_fill() +
            geom_bar(data = subdat, aes(x = month, y = monthly_aqi, fill = factor(year)), 
                     stat = "identity", position = position_dodge()) +
            labs(fill = "Year", title = paste0("Mean Monthly AQI (PM 2.5) in ", county, ", ", state)) 
            
        
    } else if (by == "week") {
        subdat <- dat %>%
            group_by(year, month, week) %>%
            summarize(weekly_aqi = mean(daily_aqi), weekly_PM2.5 = mean(daily_PM2.5))
        
        g +
            new_scale_fill() +
            geom_bar(data = subdat, aes(x = week, y = weekly_aqi, fill = factor(year)), 
                     stat = "identity", position = position_dodge()) +
            labs(fill = "Year", title = paste0("Mean Weekly AQI (PM 2.5) in", county, ", ", state))
            
        
    } else if (by == "day") {
    
        g +
            new_scale_color() +
            geom_line(data = dat, aes(x = day, y = daily_aqi, color = factor(year)), size = 1) + 
            geom_vline(xintercept = unique(dat$ed_facilities)-ymd('2020-01-01'), lty = 'dashed', color = 'darkgreen', size = 1) +
            geom_vline(xintercept = unique(dat$non_essential)-ymd('2020-01-01'), lty = 'dashed', color = 'blue', size = 1) +
            geom_vline(xintercept = unique(dat$stay_at_home)-ymd('2020-01-01'), lty = 'dashed', color = 'black', size = 1) +
            labs(color = "Year", title = paste0("Mean daily AQI in ", county, ", ", state), 
                 subtitle = paste("(green dashed line = Educational services closed,",
                                  "blue dashed line = Non-essential services closed,",
                                  "black dashed line = Stay at home order", sep = "\n"))
    }
    
}

#### Define UI ####
ui <- fluidPage(
    theme = shinytheme("flatly"),
    titlePanel("Covid-19 Projections by U.S. States and Counties"),

    sidebarLayout(
        sidebarPanel(
            tags$h4("Select state and county you are interested in:"),
            htmlOutput('state_selector'),
            htmlOutput('county_selector'),
            tags$h6("*Note: Only counties >100 cumulative cases available; models don't do well with small amount of data.")
        ), 
        
        # Output: Description, layout, and reference
        mainPanel(
            tabsetPanel(
                tabPanel("Summary plots", 
                         tags$h4('Cumulative plot'),
                         tags$div(
                             "Below are summary plots based on observed data from The New York Times. 
                             The first plot below shows cumulative cases and deaths on a", tags$b(tags$em("log")), "scale."
                             ),
                         tags$br(), 
                         plotOutput(outputId = "cum_log_plot"), 
                         tags$h4("Daily incident plots"),
                         tags$div(
                             "The plot below is for observed confirmed daily cases and deaths reported by The New York Times, but", 
                             tags$span(style="color:red", tags$b(tags$em("note that the deaths are scaled by 10"))), ", i.e., 
                             deaths are multiplied by 10 for visualization purposes.", 
                             tags$b("These are NOT actual deaths!")
                             ),
                         tags$br(),
                         plotOutput(outputId = "epi_plot")
                         ), 
                
                tabPanel("Log-linear fits",
                         tags$h4("Log-linear fit to confirmed cases data"),
                         tags$div(
                             "Below is a log-linear fit to the growth part of the incidence curve using the incidence package.
                             Log-linear models have the form", tags$br(), tags$br(), 
                             tags$em("log(y) = r * t + b"), tags$br(), tags$br(),
                             "where y is the expected incidence,", tags$br(),
                             "r is the growth rate,", tags$br(),
                             "t is the duration; and", tags$br(),
                             "b is the intercept.", 
                             tags$br(), tags$br(),
                             "Obviously you would need to fit two loglinear models -- one for the growth phase and one for the decay phase. 
                             But, here I have only fit one model to the growth phase because, in most cases for now, there are not enough 
                             data points to fit a model to the decay phase. We will be able to do that in a couple of days once we are
                             past the peak. I have written a post about this on ", 
                             tags$a(href="https://aksinghal86.github.io/covid19/posts/estimating-R-in-US/", "my blog"), "and on ",
                             tags$a(href="https://medium.com/@aksinghal.aks/exploring-covid-19-progression-in-the-us-using-r-7449cc10b0ea", "Medium."), 
                             "All the code for this analysis is available on", tags$a(href="https://github.com/aksinghal86/covid19", "my Github repo."),
                             tags$br(), tags$br()
                         ),
                         plotOutput(outputId = "loglinear_fit_plot"),
                         # Select date range to be plotted
                         fluidRow(
                             tags$em("Select start and peak dates of your choice for log-linear fit above"), tags$br(),
                             column(6, dateInput("start_date", "Start date", 
                                                 value = "2020-02-01", min = "2020-01-01", max = (today()-1))),
                             # default peak date is the one returned by the loglinear fit 
                             column(6, htmlOutput("peak_date"))
                                    
                         ),
                             "If the fit above doesn't look great, try adjusting the start date above closer to when the cases started appearing
                             and the fit may improve. It's not a super robust fit because the reported cases are low and unreliable for the 
                             most part in the US due to a lack of available testing. The model has also automatically identified a peak, 
                             but it may not be a realistic peak (again, due to underreporting of cases in the US). You can also pick a new peak that
                             you think is better!", tags$br(), tags$br(),
                             "Summary statistics are printed below. The growth rate, r, (not the same as reproductive number below)
                             should change when you changed the date above!", 
                             tags$br(), tags$br(),
                        
                        verbatimTextOutput("summary"),
                        tags$br(),
                        
                        tags$h4("Estimate of reproductive number, R0"),
                        tags$div(
                            "R0 is essentially the number of people an infected person infects. So, an R0 = 2 means that, on average, two people 
                            will catch the disease from one contagious person. We can estimate the R0 for the growth phase based on the log-linear 
                            fit above. The `lm2R0_sample` function from the `epitrix` package generates a sample of R0 values from a log-linear
                            regression like the one above."
                        ),
                        plotOutput(outputId = "r0_distr_plot"),
                        tags$div(
                            "The plot above is histogram of the sampled R0 values. Mean anticipated R0 for the county and state that you have 
                            selected is the peak of the histogram. This is pretty qualitative, but is nevertheless super helpful. You can compare
                            the R0 between different states and counties and that should give you an idea of how fast the infection is spreading
                            there relative to other places. A lower R0 may be a result of early public health interventions. For example, compare
                            Michigan and New York to California or Nevada and you'll notice that both Michigan and New York have higher R0s, 
                            which should match the news about how bad the situtation got there.",
                            tags$br(), tags$br(),
                            "As you can imagine though that R0 is more of a snapshot and doesn't really provide a good sense of reality, where 
                            you would expect the reproduction number to change over time based on public health interventions like 
                            social distancing and shelter in place. That's where the ", tags$em('effective or instantaneous 
                            reproductive number, Re,'), "comes into place.", 
                            tags$br(), tags$br(),
                        ),
                        
                        tags$h4("Estimate of effective reprodutive number, Re"),
                        tags$div(
                            "We want to measure how public health interventions in a county or a state has affected the reproductive number.
                            Using the `EpiEstim` package, we can actually estimate the instantaneous reproduction number.", 
                            tags$br(), tags$br()
                        ),
                        
                        plotOutput(outputId = "re_plot"),
                        tags$div(
                            "This is an amazingly useful plot. The bright red line marks the point at which reproduction number = 1. This 
                            is the threshold below which there would be <1 infection from each infected person meaning the epidemic has
                            been effectively brought under control and that it should die out if everything stays on course and there is 
                            no second wave.",
                            tags$br(), tags$br(), 
                            "There is another pertinent insight hiding in this plot, which is not immediately obvious. Compare the incidence
                            plot (first one) and the reproductive number (Re) peak on this plot. Notice something? While the incidence plot
                            may not indicate whether a location has reached a peak number of cases, the effective R can tell you whether the
                            number of infections have actually been going down for a while, indicating if the public health interventions
                            have been effective or not. You can read more about it here:",
                            tags$a(href = "https://medium.com/@aksinghal.aks/exploring-covid-19-progression-in-the-us-using-r-7449cc10b0ea",
                                   "this post."), tags$br(), tags$br(), 
                            "In fact, in cases where information is available on when non-essential and essential services were closed 
                            (try picking a state), you'll typically see a gradual drop-off in the effective reproductive number. Cool stuff!"
                        ),
                        tags$br(), 
                        
                        tags$h4("Force of Infection"),
                        tags$div(
                            "One last graph to consider is the 'force of infection' plot, which projects how the infections will go down if no 
                            new cases are observed.The dashed line indicates the date as of today and to the right of this line are all projections.
                            This is a confusing plot, but the thing to notice here is that if the cases are trending down", tags$em("before"), 
                            "the dashed line, then it's a good sign that the epidemic is being contained."
                        ), tags$br(),
                        plotOutput(outputId = "lambda_plot")
                ),
                
                tabPanel("Preliminary Projections", 
                         tags$br(),
                         tags$div(
                             tags$em("Preface:"), 
                             "These projections are made using logistic regression, which is pretty common in biological and
                             economic growth models because logistic regression follows a sigmoidal-shaped curve (s-curve) -- there is
                             slow growth in the beginning followed by an exponential growth phase until the resources run out, which prompts
                             a decline phase until you reach an equilibrium. Though public health interventions would slow the spread of the
                             disease, or flatten the curve so to speak, it would still be expected to follow a logistic-type growth.",
                             tags$br(), tags$br(),
                             "However, please note that these are fairly simple projections. These don't take into account any public health
                             interventions or any uncertainty or covariance within the parameters. As such, these estimates are likely to
                             be on the lower end of the spectrum. There are also no uncertainty bounds around these estimates for aforementioned
                             reasons. Nevertheless, the estimates are a really good starting point on which to improve upon.",
                             tags$br(), tags$br(),
                             "Logistic regression models follow the equation: ", tags$br(), tags$br(),
                             tags$em("C(t) = p/1 + exp(-a * (t - B))"),
                             tags$br(), tags$br(),
                             "where:", tags$br(),
                             "p is the asymptotic value of infections (the peak of the curve),", tags$br(),
                             "a is the growth rate prior to the peak of the infections, and", tags$br(),
                             "B is where rate of change is maximal.", tags$br(), tags$br(),
                             "I used the nlsLM package from the minpack.LM package. I also wrote two related posts on this on my blog",
                             tags$a(href="https://aksinghal86.github.io/covid19/posts/covid19-projections-using-logistic_regression/", "here"), 
                             "and", 
                             tags$a(href="https://aksinghal86.github.io/covid19/posts/covid19-projections-using-GRM/", "here."),
                             "All the code is also available on ", tags$a(href="https://github.com/aksinghal86/covid19", "my Github repo.")
                         ),
                         tags$br(), 

                         tags$h4("Cumulative case projections"),
                         tags$div(
                             "Even a simple logistic regression model fits the data points pretty well. The typical S-curve is also evident.
                             While the fit is pretty good, it is unable to unaccount for 'flattening the curve' measures, so it tapers out 
                             faster than it should leading to likely underestimates of actual confirmed cases."
                         ), tags$br(),
                         plotOutput(outputId = "cum_case_projections"), tags$br(),
                         
                         tags$div(
                             "According to the baseline model, below are the fitted parameters", 
                             verbatimTextOutput("cum_case_fitted_params"), 
                             "Remember, p is the total number of estimated cases (the asymptote above), a is the growth rate, and B 
                             is when the rate of change is maximal starting from the start date."
                         ), 
                         
                         tags$h4("Daily case projections"),
                         tags$div(
                             "Using the same fit above, but just uncumulated and plotted over the actual confirmed cases."
                         ),
                         plotOutput(outputId = "epi_case_projections"),
                         
                         tags$h4("Cumulative death projections"), 
                         tags$div(
                             "Now for the total number of projected deaths. Compare these to the ones made by IHME via University of Washington
                             for some context. Those are obviously vastly more robust than these estimates but these aren't too shabby either
                             for a baseline model."
                         ), tags$br(),
                         plotOutput(outputId = "cum_death_projections"), 
                         tags$div(
                             "And the fitted model parameters:",
                             verbatimTextOutput("cum_death_fitted_params")
                         ),
                         
                         tags$h4("Daily death projections"), 
                         tags$div(
                             "And now for the daily death projection based on the fit on the cumulative deaths above."
                         ), tags$br(),
                         plotOutput(outputId = "epi_death_projections"),
                         tags$br(), tags$br(), 
                         
                         tags$div(
                             "To come soon: projections building on this model that take into account various relationships between the fitted
                             parameters p, a and B!"
                         )
                     ), 
                tabPanel("AQI",
                tags$br(),
                tags$div(
                    tags$h4("AQI (work still in development)"),
                    "Below is PM2.5 AQI chart to play with for 2019 and 2020. Have not run any statistical analyses yet, 
                    so don't have any opinions to offer, except that on a perfunctory look, I don't see a whole lot of 
                    difference in air quality in most counties.", 
                    tags$br(), tags$br(), 
                    "Default values are presented by daily changes in PM2.5, but you can see changes on a weekly or 
                    monthly basis -- just select your preference in the drop down menu below. If the goal is to assess whether 
                    there has been improvement in air quality due to lockdown, it may be better to consider differences in", 
                    tags$em('monthly'), "AQI between 2019 and 2020. It would be ideal to get data for the last few years 
                    and get a feel for distribution of PM 2.5 on a daily basis rather than point estimates, but the data 
                    are rather painful to collect so I'll save that for later...",
                    # "Since most states have been in lockdown and people movement has been severely restricted, I was curious
                    # to see if there has been an improvement in air quality across the US. At first, I thought the air quality 
                    # had indeed approved, but I realized I was only looking at the 2020 data so had no context around the trends. 
                    # May be AQI always improves in the spring time due to other factors like temperature increase, etc. So, when 
                    # I started comparing AQI to previous years, lo and behold,  there was no significant difference between the 
                    # decrease in previous years and this year. This is indeed true across the board!", 
                    tags$br(), tags$br(),
                    "(Note: There are no data available for some counties in which case the graph will just show an empty space. 
                    If there are multiple sites in a county or state, those values were averaged.) ",
                    tags$br(), tags$br(),
                    fluidRow(
                        column(12, align = "center", 
                               selectInput("agg_aqi_data_by", "Aggregate data by: ", 
                                           c('day', 'week', 'month'), selected = "day")
                        )
                    ),
                    plotOutput(outputId = "aqi_plot"),
                    fluidRow(
                        column(12, align = "center", 
                               sliderInput("aqi_days", "Days since beginning of year",
                                           min = 0, max = max(aqi_data$day), value = 0)
                        )
                    ),
                    tags$br(), 
                    # "In fact the popular opinion and the narrative that the air quality has improved due to lockdowns is largely
                    # false and does not seem like an empirically derived conclusion. Even",
                    # tags$a(href="https://www.nytimes.com/2020/04/20/nyregion/coronavirus-nyc-numbers-unemployment.html",  
                    # "NY Times made the same mistake."), "Look at the last section of this article! This is actually a pretty
                    # bold example of cherry picking. They picked the poorest AQI level this year (early March) to claim a 25%
                    # improvement in AQI without taking the time to check their prior bias. Disappointing!",
                    # tags$br(), tags$br(),
                    tags$em("Preliminary observation:"), 
                    "Interestingly, traffic by itself actually does not seem to add that much to
                    improving or worsening air quality, which I definitely find hard to believe. But, the data are what they are and
                    sometimes they prove your strongest held beliefs wrong...",
                    tags$br(), tags$br(), 
                    "Data are from",
                    tags$a(href="https://www.epa.gov/outdoor-air-quality-data/download-daily-data", "U.S. EPA"),
                    "and a royal pain in the butt to fetch. I actually literally had to fetch the data for every single state
                    manually because their API sucks. Point being that these plots unfortunately won't be updated daily. The last
                    time I fetched the data was April 19, 2020; it'll probably be at least a week before this is updated to reflect
                    new data."
                    )
                ),

                tabPanel("About",
                         tags$br(),
                         tags$div(
                             "All the data used in this analysis was made available by The New York Times on their Github repo,
                             which is updated daily. Much research and inspiration from numerous blogs, websites and conversations 
                             with friends has gone behind this work.",
                             tags$br(), tags$br(),
                             "I have written some related posts on this subject matter, available on ", 
                             tags$a(href="https://aksinghal86.github.io/covid19", "my blog"), " and on ", 
                             tags$a(href="https://medium.com/@aksinghal.aks", "Medium."), "All the code is also available on ",
                             tags$a(href="https://github.com/aksinghal86/covid19", "Github."),
                             tags$br(), tags$br(), 
                             
                             tags$h4("About Ankur Singhal"),
                             "An entreprenuer by nature, a biologist by education, data scientist by career and a rock climber by passion. 
                             I have started my own consulting company, Empirical Solutions Consulting, LLC, written several scientific 
                             manuscripts, worked as a consultant for years and launched a beta version of my climbing web app", 
                             tags$em(tags$a(href="https://aksinghal86.pythonanywhere.com", "Pebblefinder.")),
                             tags$br(), tags$br(),
                             "My newest pet project is data journalism since it appears to be pretty good at combining my voracious
                             appetite for research and digging up facts, getting better at data science, and 
                             sharing with what I learn along the way with the world. I hope you will join me! :)"
                         )
                )
            )
        )
    )
)

#### Server ####
server <- function(input, output) {
    
    # Front end for sequential filtering of counties 
    output$state_selector <- renderUI({
        selectInput(inputId = 'state', 
                    label = NULL, 
                    choices = c('All states', sort(unique(state_data$state))), 
                    selected = "All states")
    })
    
    output$county_selector <- renderUI({
        counties <- county_data %>% 
            filter(state %in% input$state, cum_cases > 100) %>% 
            pull(county)
        
        counties <- c('All counties', sort(unique(counties)))
        selectInput(inputId = "county", 
                    label = NULL, 
                    choices = counties,
                    selected = 'All counties')
    })
    
    # Subset data based on state and county input by user 
    geo_loc <- reactive({
        # Filter out data
        df <- all_data %>% filter(state %in% input$state, county %in% input$county)
        
        # Error handler: an empty dataframe is passed the first time page is loaded
        # so when trying to find peak date it returns an error
        if(nrow(df) == 0) {
            peak_date <- today() + 1
        } else {
            peak_date <- create_incidence_objects(df, peak_date = today()+1)$i.peak
        }
        return(list('data' = df, 'peak_date' = peak_date))
    })

    sub_aqi_data <- reactive({
        aqi_data %>% filter(state %in% input$state, county %in% input$county, day >= input$aqi_days)
    })
    
    output$cum_log_plot <- renderPlot({
        cum_log_plot(geo_loc()$data, input$county, input$state)
    })
    
    output$epi_plot <- renderPlot({
        epi_plot(geo_loc()$data, input$county, input$state)
    })

    output$peak_date <- renderUI({
        dateInput(inputId = "peak_date", "Peak date: ", 
                  value = geo_loc()$peak_date, 
                  min = '2020-02-02', max = today())
    })
    
    # Create incidence objects for use in incidence plots (log linear fits, R0 and Re plots)
    geo_loc_by_date <- reactive({
        geo_loc()$data %>% filter(date >= input$start_date)
    })
    
    incidence_objects <- reactive({
        # Error handler
        if (is.null(input$peak_date)) {
            peak_date <- today() + 1
        } else {
            peak_date <- input$peak_date
        }
        
        create_incidence_objects(geo_loc_by_date(), peak_date)
    })
    
    output$loglinear_fit_plot <- renderPlot({
        loglinear_fit_plot(incidence_objects(), input$county, input$state)
    })
    
    output$summary <- renderPrint({
        incidence_objects()$i.fit
    })

    output$r0_distr_plot <- renderPlot({
        r0_distr_plot(incidence_objects(), input$county, input$state)
    })

    output$re_plot <- renderPlot({
        re_plot(geo_loc_by_date(), input$county, input$state)
    })

    output$lambda_plot <- renderPlot({
        lambda_plot(incidence_objects(), input$county, input$state)
    })
    
    model_fits <- reactive({
        logistic_regression_fits(geo_loc()$data, input$county, input$state)
    })
    
    output$cum_case_projections <- renderPlot({
        model_fits()$cum_cases_plot
    })
    
    output$cum_case_fitted_params <- renderPrint({
        coef(model_fits()$cases_fit)
    })

    output$epi_case_projections <- renderPlot({
        model_fits()$epi_cases_plot
    })
    
    output$cum_death_projections <- renderPlot({
        model_fits()$cum_deaths_plot
    })
    
    output$cum_death_fitted_params <- renderPrint({
        coef(model_fits()$deaths_fit)
    })
    
    output$epi_death_projections <- renderPlot({
        model_fits()$epi_deaths_plot
    })
    
    output$aqi_plot <- renderPlot({
        if(nrow(sub_aqi_data()) > 0) aqi_plot(sub_aqi_data(), input$county, input$state, by = input$agg_aqi_data_by)
    })
}

#### Run the application ####
shinyApp(ui = ui, server = server)
