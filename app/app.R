# Libraries
library(shiny)
library(shinythemes)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(incidence)
library(epitrix)
library(distcrete)
library(EpiEstim)
library(earlyR)
library(minpack.lm)



# Load data

# Public health interventions
interventions <- read.delim('data/interventions.csv', sep = ",", header = TRUE, stringsAsFactors = FALSE) %>% 
    mutate(stay_at_home = mdy(stay_at_home), 
           ed_facilities = mdy(ed_facilities), 
           non_essential = mdy(non_essential))

# COVID19 data
# state_data <- read.delim('data/covid-19-data/us-states.csv', sep=",", header = TRUE, stringsAsFactors = FALSE) %>%
state_data <- read.delim('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv',
                          sep = ',', header = TRUE, stringsAsFactors = FALSE) %>%
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


# county_data <- read.delim('data/covid-19-data/us-counties.csv', sep=",", header = TRUE, stringsAsFactors = FALSE) %>%
county_data <- read.delim('https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv',
                           sep = ',', header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(date = ymd(date), cum_cases = cases, cum_deaths = deaths) %>% 
    group_by(state, county) %>% 
    mutate(cases = c(cum_cases[1], diff(cum_cases)), 
           deaths = c(cum_deaths[1], diff(cum_deaths)))

all_data <- dplyr::bind_rows(all_us, state_data, county_data) %>% 
    left_join(interventions, by = "state") %>%
    filter(cases >= 0)

# AQI data
# aqi_data <- read.delim('../data/aqi-by-county.csv', sep = ",", header = TRUE) %>%
#     mutate(date = ymd(date_local)) %>%
#     filter(pollutant_standard == "PM25 Annual 2012") %>% 
#     select(state, state_code, county, county_code, date, aqi) 


# Define UI
ui <- fluidPage(
    theme = shinytheme("flatly"),
    titlePanel("Covid-19 Projections by U.S. States and Counties"),
                
    sidebarLayout(
        sidebarPanel(
            htmlOutput('state_selector'),
            htmlOutput('county_selector'), 
            
            # Select date range to be plotted
            dateRangeInput("date", strong("Date range"), 
                           start = "2020-02-01", end = today(), 
                           min = "2020-01-01", max = today())
        ), 
        
        # Output: Description, layout, and reference
        mainPanel(
            tabsetPanel(
                tabPanel("Summary plots", plotOutput(outputId = "cumplot"), 
                         plotOutput(outputId = "epiplot")), 
                tabPanel("Log-linear fits", 
                         tags$h4("\n\nBelow is a log-linear fit to the growth part of the incidence curve using
                         the incidence package. Log-linear models have the form log(y) = r * t + b, where y is the 
                         expected incidence, r is the growth rate and t is the duration and b is the intercept.\n
                        Try adjusting the date on the left to see how that affects the growth rate!!\n"),
                         plotOutput(outputId = "incidenceplot"), 
                         plotOutput(outputId = "r0plot")), 
                tabPanel("Effective Reproduction Rate", 
                         tags$h4("Due to a variety of factors, the incidence curves may give the perception that the
                           number of cases is growing relentlessly when the actual growth has slowed down. Looking
                           at the cumulative plot on log scale may give you a hint if it has started to flatten out.\n
                           Another way to assess this is to look at the effecive reproductive number (Re). If it's trending
                           down, it's a good sign that the number of new infections is slowing.\n
                           The effective reproductive number is based on a serial interval distribution, which I won't 
                            talk about here, but you can get more information on",
                                 tags$a(href = "https://medium.com/@aksinghal.aks/exploring-covid-19-progression-in-the-us-using-r-7449cc10b0ea",
                                        "this post.")), 
                         plotOutput(outputId = "resplot")), 
                # tabPanel("Force of infection",
                #          tags$h4("\nLastly, we can look at the force of infection plot to identify a trend. Without
                #                  worrying about the details, if the lambdas (the height of the bar), is decreasing,
                #                  it's a good sign that the number of infections is decreasing as well.\n"),
                #          plotOutput(outputId = "lambdaplot")),
                # tabPanel("AQI",
                #          plotOutput(outputId = "aqiplot"))
                tabPanel("Preliminary projections",
                         tags$h4("Projections using simple logistic regression. 
                                Note though that these are very preliminary! 
                                Haven't accounted for the uncertainty in the parameters or 
                                 accounted for public health interventions"),
                         plotOutput(outputId = "cum_case_projections"), 
                         plotOutput(outputId = "inc_case_projections"), 
                         plotOutput(outputId = "cum_death_projections"),
                         plotOutput(outputId = "inc_death_projections"))
            )
        )
    ), 
    hr(), 
    print('~~~Data made available by The New York Times~~~')
)

server <- function(input, output) {
    
    output$state_selector <- renderUI({
        selectInput(inputId = 'state', 
                    label = "Select state:", 
                    choices = c('All states', sort(unique(state_data$state))), 
                    selected = "All states")
    })
    
    output$county_selector <- renderUI({
        counties <- county_data %>% 
            filter(state %in% input$state) %>% 
            pull(county)
        
        counties <- c('All counties', sort(unique(counties)))
        selectInput(inputId = "county", 
                    label = "Select county:", 
                    choices = counties,
                    selected = 'All counties')
    })
    
    # Subset data
    selected_data <- reactive({
        all_data %>% 
            filter(state %in% input$state, 
                   county %in% input$county, 
                   date >= input$date[1], 
                   date <= input$date[2])
    })

    ########################################################################
    
    # Set default parameters for SI from Du et al. (2020)
    mu <- 3.96
    sigma <- 4.75    
    
    # Set default theme 
    my_theme <- theme_minimal() + 
        theme(legend.position = 'none',
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 18),
              plot.subtitle = element_text(size = 14),
              panel.grid.major = element_line(color = "darkgrey")
        )
    
    # Create cumulative plot
    output$cumplot <- renderPlot({
        ggplot(selected_data(), aes(x = date)) +
            geom_line(aes(y = log(cum_cases)), color = "red", size = 1) +
            geom_line(aes(y = log(cum_deaths)), color = "red", size = 1, linetype = "dashed") +
            labs(x = "", y = "N (log scale)", 
                 title = paste0("Cumulative cases and deaths in ", input$county, ", ", input$state),
                 subtitle = "(solid line = confirmed cases, dashed line = confirmed deaths)") +
            my_theme
    })
    
    # Create epi plot
    output$epiplot <- renderPlot({
        ggplot(selected_data(), aes(x = date)) +
            geom_line(aes(y = cases), color = "darkcyan",  size = 1) +
            geom_line(aes(y = deaths*10), color = "darkcyan", size = 1, linetype = "dashed") +
            annotate(geom = "text", 
                     x = input$date[1], y = max(selected_data()$cases, na.omit = TRUE), 
                     hjust = 0, vjust = 2,  
                     color = 'red', size = 5, 
                     label = 'Deaths MULITPLIED by 10\nfor visualization purposes\nNOT ACTUAL DEATHS!!!') +
            labs(x = "", y = "N", 
                 title = paste0("Confirmed daily cases in ", input$county, ", ", input$state), 
                 subtitle = "(solid line = confirmed cases, dashed line = confirmed deaths)") +
            my_theme
    })
    
    incidence_objects <- function(){
        inc_cases <- selected_data() %>%
            ungroup(state) %>%
            select(date, cases) %>%
            uncount(cases)
        
        # Create an incidence object and find the peak
        i <- incidence(inc_cases$date)
        i.peak <- incidence::find_peak(i)
        i.fit <- incidence::fit(subset(i, to = i.peak))
        
        return(list('i' = i, 
                    'i.peak' = i.peak, 
                    'i.fit' = i.fit))
    }
    
    # Create incidence plots
    output$incidenceplot <- renderPlot({
        plot(incidence_objects()$i, color="turquoise", border="white") %>%
            add_incidence_fit(incidence_objects()$i.fit) +
            labs(title = paste0("Observed and modeled incidence in ", input$county, ", ", input$state)) +
            scale_x_date(date_labels = "%b %d") +
            my_theme 
    })

    # Distribution of R0 plot
    output$r0plot <- renderPlot({
        set.seed(1)
        param <- epitrix::gamma_mucv2shapescale(mu, sigma/mu)
        w <- distcrete("gamma", interval = 1, shape = param$shape, scale = param$scale, w = 0)
        growth_R0 <- lm2R0_sample(incidence_objects()$i.fit$model, w)
        hist(growth_R0, col = "turquoise", border = "white", 
             main = "Distribution of growth R0", xlab = '') 
            my_theme 
        
    })
    
    # Instantaneous effective R plot
    output$resplot <- renderPlot({
        res <- EpiEstim::estimate_R(
            selected_data() %>%
                select(date, cases) %>%
                rename(dates = date, I = cases),
            method = "parametric_si", config = make_config(
                list(mean_si = mu, std_si = sigma)
            )
        )
        plot(res, "R", xlab = '') +
            geom_hline(yintercept = 1, color = "red", size = 1) + 
            geom_vline(xintercept = unique(selected_data()$ed_facilities), lty = 'dashed', color = 'darkgreen', size = 1) +
            geom_vline(xintercept = unique(selected_data()$non_essential), lty = 'dashed', color = 'darkblue', size = 1) + 
            geom_vline(xintercept = unique(selected_data()$stay_at_home), lty = 'dashed', color = 'darkorange', size = 1) +
            labs(x = "", title = paste("Estimated R for", input$state), 
                 subtitle = "(green dashed line = Educational services closed,\nblue dashed line = Non-essential services closed, \norange dashed line = Stay at home order)") +
            my_theme +
            theme(legend.position = 'none')
    })
    
    # Force of infection plot
    # output$lambdaplot <- renderPlot({
    #     simple_R <- get_R(incidence_objects()$i, si_mean = mu, si_sd = sigma, max_R = 8)
    #     plot(simple_R, 'lambdas', bty = 'n') +
    #         my_theme
    #     abline(v = today(), col = 'darkturquoise', lty = 'dashed', lwd = 2)
    # })
    # 
    logistic_regresson_preds <- function(){
        # set the x and y values for logistic regression
        C <- selected_data()$cum_cases
        D <- selected_data()$cum_deaths
        t <- 1:nrow(selected_data()) # nlsLM doesn't take dates
        
        # Fit logistic regression
        cases_fit <- nlsLM(C ~ p/(1 + exp(-a * (t-B))), start = list(p=10, a=0.2, B=1), lower = c(min(C), 0, 0))
        deaths_fit <- nlsLM(D ~ p/(1 + exp(-a * (t-B))), start = list(p=10, a=0.2, B=1), lower = c(min(C), 0, 0))
        
        # Predict two months in the future
        Day <- 1:(nrow(selected_data()) + 60)
        preds <- data.frame(date = selected_data()$date[1] + days(Day-1)) %>% 
            mutate(C_preds = ceiling(predict(cases_fit, list(t = Day))),
                   I_preds = c(C_preds[1], diff(C_preds)), 
                   D_preds = ceiling(predict(deaths_fit, list(t = Day))),
                   d_preds = c(D_preds[1], diff(D_preds))) %>% 
            left_join(selected_data() %>% ungroup() %>% select(date, cum_cases, cases, cum_deaths, deaths))
        
        return(preds)
    }
    
    # Logistic Regression plots
    output$cum_case_projections <- renderPlot({
        ggplot(logistic_regresson_preds(), aes(x=date)) +
            geom_line(aes(y=C_preds), color = "turquoise", size = 2) +
            geom_point(aes(y=cum_cases), color = "magenta", shape = 1, size = 3) +
            geom_line(aes(y=cum_cases), color = "gray10", linetype = "dashed") + 
            labs(x = "", y = 'Cumulative Incidence', 
                 title = paste('Cumulative cases fitted vs observed in', input$county, ",", input$state),
                 subtitle = '(turquoise = fitted incidence, magenta = observed incidence)') +
            my_theme
    })
    
    output$inc_case_projections <- renderPlot({
        ggplot(data=logistic_regresson_preds(), aes(x=date)) +
            geom_line(aes(y=I_preds), color = "turquoise", size = 2) +
            geom_point(aes(y=cases), color = "magenta", shape = 1, size = 3) +
            geom_line(aes(x=date, y=cases), color = "gray10", linetype = "dashed") + 
            labs(x = "", y = 'Daily Incidence', 
                 title = paste('Daily cases fitted vs observed in', input$county, ",", input$state),
                 subtitle = '(turquoise = fitted incidence, magenta = observed incidence)') +
            my_theme
    })
    
    
    output$cum_death_projections <- renderPlot({
        ggplot(logistic_regresson_preds(), aes(x=date)) +
            geom_line(aes(y=D_preds), color = "red", size = 2) +
            geom_point(aes(y=cum_deaths), color = "orange", shape = 1, size = 3) +
            geom_line(aes(y=cum_deaths), color = "gray10", linetype = "dashed") + 
            labs(x = "", y = 'Cumulative Incidence', 
                 title = paste('Cumulative deaths fitted vs observed in', input$county, ",", input$state),
                 subtitle = '(turquoise = fitted incidence, magenta = observed incidence)') +
            my_theme
    })
    
    output$inc_death_projections <- renderPlot({
        ggplot() +
            geom_line(data=logistic_regresson_preds(), aes(x=date, y=d_preds), color = "red", size = 2) +
            geom_point(data=selected_data(), aes(x=date, y=deaths), color = "orange", shape = 1, size = 3) +
            geom_line(data=selected_data(), aes(x=date, y=deaths), color = "gray10", linetype = "dashed") + 
            labs(x = "", y = 'Daily Incidence', 
                 title = paste('Daily deaths fitted vs observed in',  input$county, ",", input$state),
                 subtitle = '(turquoise = fitted incidence, magenta = observed incidence)') +
            my_theme
    })
    
    
    # output$aqiplot <- renderPlot({
    #     aqi_data %>% filter(state %in% input$state, county %in% input$county) %>%
    #         ggplot(., aes(x = date, y = aqi)) +
    #         geom_point(color = "turquoise") + 
    #         stat_smooth(method = "lm") +
    #         my_theme
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)
