# Do users pay attention to prices for public charging?
Data and documentation of code for scientific replication

I. Project Description

The purpose of this project is to investigate the effect of a tiered pricing policy on EV user charging behavior, using real-time data. Analysis was conducted using R (version 4.1.0), and requires the following packages:
  - Visualization and display
      - ggplot2
      - extrafont
      - reshape
  - McCrary tests
      - rddensity
      - lpdensity
  - Regression discontinuity implementation
      - rdd
      - splines
  - Quantile treatment effects
      - dplyr
      - plm
      - lmtest
      - sandwich
      - modelsummary
      - stargazer
      - lubridate
      - estimatr

All protocols were tested and validated on November 2021.

II. Data

Pseudomyzed data used for this project are contained in sessions.csv, which encompasses 22052 charging sessions at Georgia Tech EV parking lots from July 5, 2014 to December 31, 2019. The data were downloaded using the ChargePoint API, with the API key and password provided to the authors by Georgia Tech Parking and Transportation Services. Each observation represents a charging session. Only valid sessions with a start and end time were recorded.

The available list of fields in the derived data to support replication includes:
  - Start Time, End Time of session (including breakdowns by year, month, day, hour, minute, and second)
  - Energy Consumed (kWh per session)
  - Station ID and Port Number
  - Pseudomyzed Session ID
  - Pseudomyzed User ID
  - Total duration of session
  - Parking Lot

No personally identifiable information was collected or shared.

For further details and instructions, see the comments in the R file EV_charging_project.R.















