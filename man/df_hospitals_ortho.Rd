\name{df_hospitals_ortho}
\alias{df_hospitals_ortho}
\docType{data}
\title{
Orthopedic Device Sales
}
\description{
Sales data of orthopedic device company to client hospitals over almost 2 years. 15 variables, 4703 hospitals.  Unit of observation is a unique hospital.
}
\usage{data(df_hospitals_ortho)}
\format{
  Variables are:
  
\code{zip}:               'character': Postal code.

\code{hid}:               'character': Hospital ID.

\code{city}:              'character': Hospital city.

\code{state}:             'character': Hospital state.

\code{tot_sales}:         'numeric': Total sales to hospital.

\code{tot_knee}:          'numeric': Number of knee operations.

\code{tot_hip}:           'numeric': Number of hip operations.

\code{beds}:              'numeric': Number of beds.

\code{rehab_beds}:        'numeric': Number of beds dedicated for rehabilitation.

\code{outpatient_visits}: 'numeric': Number of outpatient visits.

\code{adm_costs}:         'numeric': Administrative costs ($1000's / yr).

\code{revenue_inpatient}: 'numeric': Inpatient revenue.

\code{is_teaching}:       'numeric': Is teaching hospital?

\code{has_trauma}:        'numeric': Has trauma center?

\code{has_rehab}:         'numeric': Offers rehabilitation?
 }

\details{
  This data frame has attribute `modelvars' which gives names of numeric model variables.
}

\source{
Data adapted from `c84.dat' from Statistical Consulting, Javier Cabrera and Andrew McDougall.
}
\references{
Statistical Consulting, Javier Cabrera and Andrew McDougall.  Springer, Piscataway, NJ, 2002.
}
\examples{
data(df_hospitals_ortho)

tail(df_hospitals_ortho)

dim(df_hospitals_ortho)

attr(df_hospitals_ortho, "modelvars")
}
\keyword{datasets}
