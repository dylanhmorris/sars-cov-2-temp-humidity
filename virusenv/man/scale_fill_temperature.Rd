% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_style.R
\name{scale_fill_temperature}
\alias{scale_fill_temperature}
\title{scale_fill_temperature}
\usage{
scale_fill_temperature(
  limits = c(0, 40),
  end = 0.9,
  name = "temperature (°C)",
  ...
)
}
\arguments{
\item{limits}{limits of the temperature scale (default 0 and 40)}

\item{end}{endpoint of the colorscale within the viridis palette
(default 0.9, to avoid invisible hot temperatures)}

\item{name}{display name for the scale, default temperature
(deg. C)}

\item{...}{additional parameters passed to scale_fill_viridis()}
}
\value{
fill scale
}
\description{
default fill scale for temperature
for this project: uses the viridis
magma palette
}
