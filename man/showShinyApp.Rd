\name{showShinyApp}
\alias{showShinyApp}
\title{
Display a 'shiny' application
}
\description{
Displays one of the built-in interactive 'shiny' applications in the browser. See Details for the apps available.
}
\usage{
showShinyApp(topic)
}

\arguments{
  \item{topic}{
The name of the shiny app to display.
}

}
\details{
Two apps are currently included in the \pkg{wiqid} package:

\emph{"Beta"} displays a beta distribution and sliders which allow you to change the parameters. You can also input binomial data and obtain the conjugate beta posterior distribution.

\emph{"Gamma"} displays a gamma distribution with variable parameters, and can produce the conjugate gamma posterior for Poisson-distributed count data.
}
\value{
Nothing useful. The function is run for its side effect.
}

\author{
A much simplified version of code by Jason Bryer on Github at https://github.com/jbryer/IS606, adapted by Mike Meredith. 
}

\examples{
showShinyApp() # Shows a list of available apps

# The following command will cause R to freeze until you press Esc:
# showShinyApp("Beta")

}